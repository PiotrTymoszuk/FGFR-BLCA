# Fetches images with IHC staining from the Human Protein Atlas

  insert_head()

# downloading -------

  insert_msg('Downloading')

  hpa$data[, c("image_url", "image_path")] %>%
    set_names(c('url', 'file')) %>%
    pwalk(download_image)

# END -------

  insert_tail()
